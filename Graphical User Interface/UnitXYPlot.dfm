object FormXYPlot: TFormXYPlot
  Left = 0
  Top = 0
  Caption = 'XYPlot parameters'
  ClientHeight = 188
  ClientWidth = 277
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  OnCreate = FormCreate
  PixelsPerInch = 96
  TextHeight = 13
  object GroupBoxXYPlot: TGroupBox
    Left = 8
    Top = 8
    Width = 265
    Height = 177
    Caption = 'XY Plot'
    TabOrder = 0
    object LabelXo: TLabel
      Left = 22
      Top = 27
      Width = 12
      Height = 13
      Caption = 'Xo'
    end
    object LabelYo: TLabel
      Left = 22
      Top = 54
      Width = 12
      Height = 13
      Caption = 'Yo'
    end
    object LabelZo: TLabel
      Left = 22
      Top = 80
      Width = 12
      Height = 13
      Caption = 'Zo'
    end
    object Labeldirectional: TLabel
      Left = 8
      Top = 120
      Width = 49
      Height = 13
      Caption = 'directional'
    end
    object LabelXdim: TLabel
      Left = 176
      Top = 32
      Width = 47
      Height = 13
      Caption = 'LabelXdim'
    end
    object LabelYdim: TLabel
      Left = 176
      Top = 56
      Width = 47
      Height = 13
      Caption = 'LabelYdim'
    end
    object LabelZdim: TLabel
      Left = 176
      Top = 80
      Width = 47
      Height = 13
      Caption = 'LabelZdim'
    end
    object ButtonApply: TButton
      Left = 128
      Top = 144
      Width = 75
      Height = 25
      Caption = 'Apply'
      TabOrder = 0
      OnClick = ButtonApplyClick
    end
    object EditXo: TEdit
      Left = 40
      Top = 24
      Width = 121
      Height = 21
      TabOrder = 1
    end
    object EditYo: TEdit
      Left = 40
      Top = 51
      Width = 121
      Height = 21
      TabOrder = 2
    end
    object EditZo: TEdit
      Left = 40
      Top = 78
      Width = 121
      Height = 21
      TabOrder = 3
    end
    object ComboBoxdirectional: TComboBox
      Left = 63
      Top = 117
      Width = 145
      Height = 21
      ItemIndex = 0
      TabOrder = 4
      Text = 'X'
      Items.Strings = (
        'X'
        'Y'
        'Z')
    end
  end
end
