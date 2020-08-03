object Formamg1r5Parameters: TFormamg1r5Parameters
  Left = 0
  Top = 0
  AutoSize = True
  Caption = 'amg1r5 Parameters'
  ClientHeight = 237
  ClientWidth = 273
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object Labelstrongthreshold: TLabel
    Left = 16
    Top = 176
    Width = 80
    Height = 13
    Caption = 'Strong threshold'
  end
  object LabelF2F: TLabel
    Left = 152
    Top = 176
    Width = 28
    Height = 13
    Caption = 'F to F'
  end
  object GroupBoxStabilisation: TGroupBox
    Left = 0
    Top = 0
    Width = 273
    Height = 57
    Caption = 'Stabilisation'
    TabOrder = 0
    object ComboBoxStabilization: TComboBox
      Left = 3
      Top = 24
      Width = 198
      Height = 21
      ItemIndex = 1
      TabOrder = 0
      Text = 'BiCGStab'
      Items.Strings = (
        'none'
        'BiCGStab'
        'FGMRes'
        'NonLinear')
    end
  end
  object ButtonApply: TButton
    Left = 186
    Top = 212
    Width = 75
    Height = 25
    Caption = 'Apply'
    TabOrder = 1
    OnClick = ButtonApplyClick
  end
  object CheckBox_amg1r6: TCheckBox
    Left = 0
    Top = 212
    Width = 59
    Height = 17
    Caption = 'amg1r6'
    TabOrder = 2
  end
  object GroupBox1: TGroupBox
    Left = 0
    Top = 67
    Width = 261
    Height = 86
    Caption = 'number of smoothers steps'
    TabOrder = 3
    object Labelpre: TLabel
      Left = 3
      Top = 24
      Width = 51
      Height = 13
      Caption = 'presmooth'
    end
    object LabelpostSmooth: TLabel
      Left = 3
      Top = 56
      Width = 56
      Height = 13
      Caption = 'postsmooth'
    end
    object ComboBoxNumber_of_smootherssteps: TComboBox
      Left = 69
      Top = 24
      Width = 29
      Height = 21
      ItemIndex = 0
      TabOrder = 0
      Text = '1'
      Items.Strings = (
        '1'
        '2'
        '3')
    end
    object ComboBoxTypeSmoother: TComboBox
      Left = 104
      Top = 24
      Width = 145
      Height = 21
      ItemIndex = 0
      TabOrder = 1
      Text = 'C/ F relaxation'
      Items.Strings = (
        'C/ F relaxation'
        'FULL GS SWEEP'
        'FULL MORE COLOR SWEEP'
        'ilu (lfil) Incompleate LU')
    end
    object ComboBoxNumber_of_post_smooth: TComboBox
      Left = 70
      Top = 51
      Width = 28
      Height = 21
      ItemIndex = 0
      TabOrder = 2
      Text = '1'
      Items.Strings = (
        '1'
        '2'
        '3')
    end
    object ComboBoxTypePostSmoother: TComboBox
      Left = 104
      Top = 51
      Width = 145
      Height = 21
      ItemIndex = 0
      TabOrder = 3
      Text = 'C/ F relaxation'
      Items.Strings = (
        'C/ F relaxation'
        'FULL GS SWEEP'
        'FULL MORE COLOR SWEEP'
        'ilu (lfil) Incompleate LU')
    end
  end
  object RGthresholds: TRadioGroup
    Left = 0
    Top = 159
    Width = 261
    Height = 47
    Caption = 'thresholds'
    TabOrder = 4
  end
  object Editstrongthreshold: TEdit
    Left = 102
    Top = 173
    Width = 43
    Height = 21
    TabOrder = 5
    Text = '0.25'
  end
  object EditF2F: TEdit
    Left = 186
    Top = 173
    Width = 47
    Height = 21
    TabOrder = 6
    Text = '0.35'
  end
end
