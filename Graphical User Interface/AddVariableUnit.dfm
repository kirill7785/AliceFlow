object AddVariableForm: TAddVariableForm
  Left = 354
  Top = 150
  BorderIcons = []
  Caption = 'Add Variable'
  ClientHeight = 166
  ClientWidth = 269
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object pnlbase: TPanel
    Left = 8
    Top = 8
    Width = 257
    Height = 153
    Color = clMoneyGreen
    TabOrder = 0
    object lblnamevar: TLabel
      Left = 16
      Top = 24
      Width = 27
      Height = 13
      Caption = 'Name'
    end
    object lblname: TLabel
      Left = 56
      Top = 24
      Width = 3
      Height = 13
    end
    object lblvalue: TLabel
      Left = 16
      Top = 64
      Width = 26
      Height = 13
      Caption = 'Value'
    end
    object edtvalue: TEdit
      Left = 56
      Top = 56
      Width = 185
      Height = 21
      TabOrder = 0
    end
    object btnApply: TButton
      Left = 136
      Top = 104
      Width = 75
      Height = 25
      Caption = 'Apply'
      TabOrder = 1
      OnClick = btnApplyClick
    end
  end
end
