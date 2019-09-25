object FormAMGCLParameters: TFormAMGCLParameters
  Left = 0
  Top = 0
  Caption = 'AMGCL'
  ClientHeight = 325
  ClientWidth = 194
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object Label1: TLabel
    Left = 8
    Top = 296
    Width = 90
    Height = 13
    Caption = 'Denis Demidov  alg'
  end
  object RadioGroupAMGCLsmoother1: TRadioGroup
    Left = 0
    Top = 0
    Width = 185
    Height = 145
    Caption = 'AMGCL smoother'
    Color = clMoneyGreen
    Columns = 2
    ItemIndex = 0
    Items.Strings = (
      'spai0'
      'ilu0'
      'gauss-seidel'
      'damped-jacobi'
      'spai1'
      'chebyshev'
      'iluk, k=1'
      'iluk, k=2')
    ParentBackground = False
    ParentColor = False
    TabOrder = 0
  end
  object Button1: TButton
    Left = 109
    Top = 291
    Width = 75
    Height = 25
    Caption = 'Apply'
    TabOrder = 1
    OnClick = Button1Click
  end
  object RadioGroupAMGCLCoarseningType: TRadioGroup
    Left = -1
    Top = 151
    Width = 185
    Height = 99
    Caption = 'AMGCL coarsening type'
    ItemIndex = 1
    Items.Strings = (
      'Ruge-Stueben amg1r5 analog'
      'smoother aggregation')
    TabOrder = 2
  end
  object ComboBoxIterator: TComboBox
    Left = 8
    Top = 264
    Width = 169
    Height = 21
    ItemIndex = 0
    TabOrder = 3
    Text = 'BiCGStab'
    Items.Strings = (
      'BiCGStab'
      'FGMRES')
  end
end
